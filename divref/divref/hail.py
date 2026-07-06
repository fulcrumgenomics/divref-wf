import os
from pathlib import Path

import hail as hl
import pyspark


def _export_gcs_credentials(gcs_credentials_path: Path | None) -> None:
    """
    Export GOOGLE_APPLICATION_CREDENTIALS for the JVM subprocess (GCS mode only).

    Args:
        gcs_credentials_path: Path to the ADC JSON file; must be provided and must exist.

    Raises:
        ValueError: If `gcs_credentials_path` is None.
        FileNotFoundError: If the file does not exist.
    """
    if gcs_credentials_path is None:
        raise ValueError("gcs_credentials_path is required when use_s3 is False.")
    if not gcs_credentials_path.is_file():
        raise FileNotFoundError(
            f"GCS credentials file not found at {gcs_credentials_path}. Run "
            "`gcloud auth application-default login` or pass a valid --gcs-credentials-path."
        )
    if "GOOGLE_APPLICATION_CREDENTIALS" not in os.environ:
        os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = str(gcs_credentials_path)


def hail_init(
    *,
    gcs_credentials_path: Path | None = None,
    spark_driver_memory_gb: int = 1,
    spark_executor_memory_gb: int = 1,
    use_s3: bool = False,
) -> None:
    """
    Initialize Hail with either the GCS connector or the S3A connector.

    When `use_s3` is `False` (the default), sets `GOOGLE_APPLICATION_CREDENTIALS`
    so the JVM subprocess inherits it, then starts Hail with the GCS connector JAR on
    the Spark classpath. When `use_s3` is `True`, the S3A connector JARs
    (`hadoop-aws` and `aws-java-sdk-bundle`) are loaded and the S3A Spark configs
    are set instead — the GCS connector is not required. S3 reads use
    `AnonymousAWSCredentialsProvider` because every input the workflow consumes
    (`gnomad-public-us-east-1`, `broad-references`) is on a public Open Data
    bucket that allows anonymous reads. No AWS credentials are needed.

    Args:
        gcs_credentials_path: Absolute path to a GCP Application Default Credentials
            JSON file. Required (and must exist) when `use_s3` is `False`; ignored
            otherwise. When `GOOGLE_APPLICATION_CREDENTIALS` is not already set, it is
            exported to the environment before Hail starts.
        spark_driver_memory_gb: Memory in GB to allocate to the Spark driver.
        spark_executor_memory_gb: Memory in GB to allocate to the Spark executor.
        use_s3: If `True`, validate and load the S3A connector JARs and configure
            S3A Spark properties (the GCS connector is skipped). Leave `False` for
            GCS workloads.

    Raises:
        ValueError: If `spark_driver_memory_gb` or `spark_executor_memory_gb`
            is less than 1, or if `use_s3` is `False` and `gcs_credentials_path`
            is `None`.
        FileNotFoundError: If `use_s3` is `False` and either the credentials file at
            `gcs_credentials_path` or the GCS connector JAR is missing, or if `use_s3`
            is `True` and either S3A connector JAR is missing.
    """
    if spark_driver_memory_gb < 1:
        raise ValueError(
            f"Spark driver memory must be at least 1GB. Saw {spark_driver_memory_gb}GB."
        )
    if spark_executor_memory_gb < 1:
        raise ValueError(
            f"Spark executor memory must be at least 1GB. Saw {spark_executor_memory_gb}GB."
        )

    os.environ["PYSPARK_SUBMIT_ARGS"] = (
        f"--driver-memory {spark_driver_memory_gb}g "
        f"--executor-memory {spark_executor_memory_gb}g "
        "pyspark-shell"
    )

    if not use_s3:
        _export_gcs_credentials(gcs_credentials_path)

    jars_dir = Path(pyspark.__path__[0]) / "jars"
    cloud_jars: list[str] = []
    spark_conf: dict[str, str] = {}

    if use_s3:
        hadoop_aws_jar = jars_dir / "hadoop-aws.jar"
        aws_sdk_bundle_jar = jars_dir / "aws-java-sdk-bundle.jar"
        for jar in (hadoop_aws_jar, aws_sdk_bundle_jar):
            if not jar.exists():
                raise FileNotFoundError(
                    f"S3 connector JAR not found at {jar}. Run 'pixi run setup-s3' to download it."
                )
        cloud_jars.extend([str(hadoop_aws_jar), str(aws_sdk_bundle_jar)])
        spark_conf.update({
            "spark.hadoop.fs.s3a.impl": "org.apache.hadoop.fs.s3a.S3AFileSystem",
            # All workflow inputs live on public Open Data buckets that allow anonymous
            # reads; this avoids needing credentials and bypasses any restrictive IAM role
            # that might be present on the host.
            "spark.hadoop.fs.s3a.aws.credentials.provider": "org.apache.hadoop.fs.s3a.AnonymousAWSCredentialsProvider",  # noqa: E501
            # Random-read optimizations for Hail/Parquet workloads. S3A defaults to
            # sequential fadvise, which re-opens the HTTP connection on every backward
            # seek; Hail does many seeks per partition, so `random` is much faster.
            "spark.hadoop.fs.s3a.experimental.input.fadvise": "random",
            "spark.hadoop.fs.s3a.readahead.range": "64K",
            # Defaults (15 / 10) bottleneck partition-parallel reads.
            "spark.hadoop.fs.s3a.connection.maximum": "200",
            "spark.hadoop.fs.s3a.threads.max": "64",
        })
    else:
        gcs_jar = jars_dir / "gcs-connector.jar"
        if not gcs_jar.is_file():
            raise FileNotFoundError(
                f"GCS connector JAR not found at {gcs_jar}. "
                "Run 'pixi run setup-gcs' to download it."
            )
        cloud_jars.append(str(gcs_jar))
        spark_conf.update({
            "spark.hadoop.fs.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem",
            "spark.hadoop.fs.AbstractFileSystem.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS",  # noqa: E501
        })

    cloud_jars_str = ",".join(cloud_jars)
    spark_conf["spark.jars"] = cloud_jars_str
    spark_conf["spark.driver.extraClassPath"] = cloud_jars_str

    hl.init(spark_conf=spark_conf)
