#!/usr/bin/env python3

from pathlib import Path

from snakemake.api import (
    OutputSettings,
    ResourceSettings,
    WorkflowSettings,
    SnakemakeApi,
    StorageSettings,
    ConfigSettings,
    DAGSettings,
    ExecutionSettings
)


def run_snakemake(snake_config,my_snakefile,v,config):
    pshell = False
    if v:
        pshell = True
    with SnakemakeApi(
        OutputSettings(
            printshellcmds=pshell,
            verbose=v
            # log_handler_settings=custom_logger,
        )
    ) as snakemake_api:
        try:
        
            workflow_api = snakemake_api.workflow(
                resource_settings=ResourceSettings(
                    cores=config[KEY_THREADS]
                    ),
                config_settings=ConfigSettings(
                    config=config
                ),
                snakefile=Path(my_snakefile),
                workdir=Path(config[KEY_TEMPDIR]),
                )
            dag_api = workflow_api.dag(
                dag_settings=DAGSettings(
                    forceall=True,
                    force_incomplete=True
                ),
                    )
            dag_api.execute_workflow(
                execution_settings=ExecutionSettings(
                    lock=False,
                    keep_incomplete=True
                ),
            )


        except Exception as e:
            snakemake_api.print_exception(e)
            return False
    
    return True

