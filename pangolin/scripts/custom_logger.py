
""" 
Adapted from snakemake source code
__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015-2019, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"
"""

import logging as _logging
import platform
import time
import datetime
import sys
import os
import json
import threading
import tempfile
from functools import partial
import inspect
import textwrap

from snakemake.common import DYNAMIC_FILL
from snakemake.common import Mode


class ColorizingStreamHandler(_logging.StreamHandler):

    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
    RESET_SEQ = "\033[0m"
    COLOR_SEQ = "\033[%dm"
    BOLD_SEQ = "\033[1m"

    colors = {
        "WARNING": YELLOW,
        "INFO": GREEN,
        "DEBUG": BLUE,
        "CRITICAL": RED,
        "ERROR": RED,
    }

    def __init__(
        self, nocolor=False, stream=sys.stderr, use_threads=False, mode=Mode.default
    ):
        super().__init__(stream=stream)

        self._output_lock = threading.Lock()

        self.nocolor = nocolor or not self.can_color_tty(mode)

    def can_color_tty(self, mode):
        if "TERM" in os.environ and os.environ["TERM"] == "dumb":
            return False
        if mode == Mode.subprocess:
            return True
        return self.is_tty and not platform.system() == "Windows"

    @property
    def is_tty(self):
        isatty = getattr(self.stream, "isatty", None)
        return isatty and isatty()

    def emit(self, record):
        with self._output_lock:
            try:
                self.format(record)  # add the message to the record
                self.stream.write(self.decorate(record))
                self.stream.write(getattr(self, "terminator", "\n"))
                self.flush()
            except BrokenPipeError as e:
                pass
            except (KeyboardInterrupt, SystemExit):
                # ignore any exceptions in these cases as any relevant messages have been printed before
                pass
            except Exception as e:
                pass

    def decorate(self, record):
        message = record.message
        message = [message]
        if not self.nocolor and record.levelname in self.colors:
            message.insert(0, self.COLOR_SEQ % (30 + self.colors[record.levelname]))
            message.append(self.RESET_SEQ)
        return "".join(message)


class SlackLogger:
    def __init__(self):
        from slacker import Slacker

        self.token = os.getenv("SLACK_TOKEN")
        if not self.token:
            print(
                "The use of slack logging requires the user to set a user specific slack legacy token to the SLACK_TOKEN environment variable. Set this variable by 'export SLACK_TOKEN=your_token'. To generate your token please visit https://api.slack.com/custom-integrations/legacy-tokens."
            )
            exit(-1)
        self.slack = Slacker(self.token)
        # Check for success
        try:
            auth = self.slack.auth.test().body
        except Exception:
            print(
                "Slack connection failed. Please compare your provided slack token exported in the SLACK_TOKEN environment variable with your online token at https://api.slack.com/custom-integrations/legacy-tokens. A different token can be set up by 'export SLACK_TOKEN=your_token'."
            )
            exit(-1)
        self.own_id = auth["user_id"]
        self.error_occured = False

    def log_handler(self, msg):
        if msg["level"] == "error" and not self.error_occured:
            self.slack.chat.post_message(
                self.own_id, text="At least one error occured.", username="snakemake"
            )
            self.error_occured = True

        if msg["level"] == "progress" and msg["done"] == msg["total"]:
            # workflow finished
            self.slack.chat.post_message(
                self.own_id, text="Workflow complete.", username="snakemake"
            )


class Logger:
    def __init__(self):
        self.logger = _logging.getLogger(__name__)
        self.log_handler = [self.text_handler]
        self.stream_handler = None
        self.printshellcmds = False
        self.printreason = False
        self.debug_dag = False
        self.quiet = False
        self.logfile = None
        self.last_msg_was_job_info = False
        self.mode = Mode.default
        self.show_failed_logs = False
        self.logfile_handler = None

    def setup_logfile(self):
        if self.mode == Mode.default:
            os.makedirs(os.path.join(".snakemake", "log"), exist_ok=True)
            self.logfile = os.path.abspath(
                os.path.join(
                    ".snakemake",
                    "log",
                    datetime.datetime.now().isoformat().replace(":", "")
                    + ".snakemake.log",
                )
            )

            self.logfile_handler = _logging.FileHandler(self.logfile)
            self.logger.addHandler(self.logfile_handler)

    def cleanup(self):
        if self.mode == Mode.default and self.logfile_handler is not None:
            self.logger.removeHandler(self.logfile_handler)
            self.logfile_handler.close()
        self.log_handler = [self.text_handler]

    def get_logfile(self):
        if self.logfile is not None:
            self.logfile_handler.flush()
        return self.logfile

    def remove_logfile(self):
        if self.mode == Mode.default:
            self.logfile_handler.close()
            os.remove(self.logfile)

    def handler(self, msg):
        for handler in self.log_handler:
            handler(msg)

    def set_stream_handler(self, stream_handler):
        if self.stream_handler is not None:
            self.logger.removeHandler(self.stream_handler)
        self.stream_handler = stream_handler
        self.logger.addHandler(stream_handler)

    def set_level(self, level):
        self.logger.setLevel(level)

    def logfile_hint(self):
        if self.mode == Mode.default:
            logfile = self.get_logfile()
            self.info("Complete log: {}".format(logfile))

    def location(self, msg):
        callerframerecord = inspect.stack()[1]
        frame = callerframerecord[0]
        info = inspect.getframeinfo(frame)
        self.debug(
            "{}: {info.filename}, {info.function}, {info.lineno}".format(msg, info=info)
        )

    def info(self, msg, indent=False):
        self.handler(dict(level="info", msg=msg, indent=indent))

    def warning(self, msg):
        self.handler(dict(level="warning", msg=msg))

    def debug(self, msg):
        self.handler(dict(level="debug", msg=msg))

    def error(self, msg):
        self.handler(dict(level="error", msg=msg))

    def progress(self, done=None, total=None):
        self.handler(dict(level="progress", done=done, total=total))

    def resources_info(self, msg):
        self.handler(dict(level="resources_info", msg=msg))

    def run_info(self, msg):
        self.handler(dict(level="run_info", msg=msg))

    def group_info(self, **msg):
        msg["level"] = "group_info"
        self.handler(msg)

    def job_info(self, **msg):
        msg["level"] = "job_info"
        self.handler(msg)

    def job_error(self, **msg):
        msg["level"] = "job_error"
        self.handler(msg)

    def group_error(self, **msg):
        msg["level"] = "group_error"
        self.handler(msg)

    def dag_debug(self, msg):
        self.handler(dict(level="dag_debug", **msg))

    def shellcmd(self, msg, indent=False):
        if msg is not None:
            msg = dict(level="shellcmd", msg=msg)
            msg["indent"] = indent
            self.handler(msg)

    def job_finished(self, **msg):
        msg["level"] = "job_finished"
        self.handler(msg)

    def rule_info(self, **msg):
        msg["level"] = "rule_info"
        self.handler(msg)

    def d3dag(self, **msg):
        msg["level"] = "d3dag"
        self.handler(msg)

    def text_handler(self, msg):
        """The default snakemake log handler.
        Prints the output to the console.
        Args:
            msg (dict):     the log message dictionary
        """
        
        def job_info(msg):
            def format_item(item, omit=None, valueformat=str):
                value = msg[item]
                if value != omit:
                    return "".format(item, valueformat(value))

            yield "{}{} {}:".format(
                "local" if msg["local"] else "",
                "checkpoint" if msg["is_checkpoint"] else "rule",
                msg["name"],
            )
            for item in ["input", "output", "log"]:
                fmt = format_item(item, omit=[], valueformat=", ".join)
                if fmt != None:
                    yield fmt

            singleitems = ["jobid", "benchmark"]
            if self.printreason:
                singleitems.append("reason")
            for item in singleitems:
                fmt = format_item(item, omit=None)
                if fmt != None:
                    yield fmt

            wildcards = format_wildcards(msg["wildcards"])

            for item, omit in zip("priority threads".split(), [0, 1]):
                fmt = format_item(item, omit=omit)
                if fmt != None:
                    yield fmt

            resources = format_resources(msg["resources"])

        def indent(item):
            if msg.get("indent"):
                return "    " + item
            else:
                return item

        def timestamp():
            self.logger.info(indent("[{}]".format(time.asctime())))
        
        level = msg["level"]
        if level == "job_error":
            timestamp()
            # self.logger.error(indent("Error in {}:".format(msg["name"])))
            
            # if msg["log"]:
            #     self.logger.error(
            #         indent(
            #             "    log: {} (check log file(s) for error message)".format(
            #                 ", ".join(msg["log"])
            #             )
            #         )
            #     )
            # if msg["conda_env"]:
            #     self.logger.error(indent("    conda-env: {}".format(msg["conda_env"])))
            # if msg["shellcmd"]:
            #     self.logger.error(
            #         indent(
            #             "Error: the following command exited with an error\n{}".format(
            #                 msg["shellcmd"]
            #             )
            #         )
            #     )

            for item in msg["aux"].items():
                self.logger.error(indent("    {}: {}".format(*item)))

            if self.show_failed_logs and msg["log"]:
                for f in msg["log"]:
                    try:
                        self.logger.error("Logfile {}:\n{}".format(f, open(f).read()))
                    except FileNotFoundError:
                        self.logger.error("Logfile {} not found.".format(f))

            self.logger.error("")
        elif level == "group_error":
            timestamp()
            self.logger.error("Error in group job {}:".format(msg["groupid"]))
        else:
            if level == "error":
                # self.logger.error(msg["msg"])
                pass
            elif level == "debug":
                self.logger.debug(msg["msg"])
            elif level == "job_finished" and not self.quiet:
                timestamp()
                self.logger.info("Finished job {}.".format(msg["jobid"]))
                pass
            elif level == "d3dag":
                print(json.dumps({"nodes": msg["nodes"], "links": msg["edges"]}))
            elif level == "dag_debug":
                if self.debug_dag:
                    if "file" in msg:
                        # self.logger.warning(
                        #     "file {file}:\n    {msg}\n{exception}".format(
                        #         file=msg["file"],
                        #         msg=msg["msg"],
                        #         exception=textwrap.indent(
                        #             str(msg["exception"]), "    "
                        #         ),
                        #     )
                        # )
                        pass
                    else:
                        job = msg["job"]

            self.last_msg_was_job_info = True


def format_dict(dict_like, omit_keys=[], omit_values=[]):
    from snakemake.io import Namedlist

    if isinstance(dict_like, Namedlist):
        items = dict_like.items()
    elif isinstance(dict_like, dict):
        items = dict_like.items()
    else:
        raise ValueError(
            "bug: format_dict applied to something neither a dict nor a Namedlist"
        )
    return ", ".join(
        "{}={}".format(name, str(value))
        for name, value in items
        if name not in omit_keys and value not in omit_values
    )


format_resources = partial(format_dict, omit_keys={"_cores", "_nodes"})
format_wildcards = partial(format_dict, omit_values={DYNAMIC_FILL})


def format_resource_names(resources, omit_resources="_cores _nodes".split()):
    return ", ".join(name for name in resources if name not in omit_resources)


logger = Logger()


def setup_logger(
    handler=[],
    quiet=True,
    printshellcmds=False,
    printreason=False,
    debug_dag=False,
    nocolor=False,
    stdout=False,
    debug=False,
    use_threads=False,
    mode=Mode.default,
    show_failed_logs=False,
):
    logger.log_handler.extend(handler)

    # console output only if no custom logger was specified
    stream_handler = ColorizingStreamHandler(
        nocolor=nocolor,
        stream=sys.stdout if stdout else sys.stderr,
        use_threads=use_threads,
        mode=mode,
    )
    logger.set_stream_handler(stream_handler)

    logger.set_level(_logging.DEBUG if debug else _logging.INFO)
    logger.quiet = True
    logger.printshellcmds = False
    logger.printreason = printreason
    logger.debug_dag = debug_dag
    logger.mode = mode
    logger.show_failed_logs = show_failed_logs
