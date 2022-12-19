"""module for logging construction of the SFR package
Based on the logger.py module in the pyemu project
Copyright (c) 2014
https://github.com/jtwhite79/pyemu
"""
from __future__ import print_function, division
import os
from pathlib import Path
import time
from datetime import datetime
import io
from contextlib import redirect_stdout
import warnings
import copy
import sfrmaker


class Logger(object):
    """ A basic class for logging events
        if filename is passed, then a file handle is opened.

    Parameters
    ----------
    filename : str or file handle
        File to write logged events to.
        by default, 'sfrmaker.logger'
    mode : str {'w', 'a'}
        Option to make a new file ('w') or append to an
        existing logger file ('a').
    echo : bool
        Flag to cause logged events to be echoed to the screen.

    Attributes
    ----------
    filename : str
        Path to the logger file
    f : file handle
        File handle that the logger being written to
    items : dict
        Dictionary holding events to be logged.  If a logger entry is
        not in `items`, then it is treated as a new entry with the string
        being the key and the datetime as the value.  If a logger entry is
        in `items`, then the end time and delta time are written and
        the item is popped from the keys.

    """
    def __init__(self, filename='sfrmaker.logger', mode='w', echo=False):
        defaults = {''}
        self.items = {}
        self.echo = bool(echo)
        if isinstance(filename, str) or isinstance(filename, Path):
            self.filename = filename
            if os.path.exists(filename) and mode == 'a':
                self.f = open(filename, 'a')
            else:
                self.f = open(filename, 'w')
                self.statement('SFRmaker v. {}'.format(sfrmaker.__version__),
                               log_time=False)
        # work with an open file handle
        else:
            self.f = filename
            self.filename = filename.name
            self.mode = filename.mode
        self.t = datetime.now()

    def statement(self, phrase,
                  log_time=True, echo=None):
        """ logger a one time statement

        Parameters
        ----------
        phrase : str
            statement to logger

        """
        t = ''
        if log_time:
            t = str(datetime.now()) + ' '
        s = t + str(phrase) + '\n'
        echo = self.echo if echo is None else echo
        if echo:
            print(s, end='')
        if self.filename:
            self.f.write(s)
            self.f.flush()

    def log_package_version(self, pkg):
        if isinstance(pkg, str):
            pkg = __import__(pkg)
        self.statement('{} v {}'.format(pkg.__name__,
                                        pkg.__version__,
                                        log_time=False))

    def log_file_and_date_modified(self, filepath, prefix='', suffix=''):
        txt = '{}, last modified: {}'.format(filepath,
                                             time.ctime(os.path.getmtime(filepath)))
        self.statement(prefix + txt + suffix,
                       log_time=False)

    def log_fn_w_stdout(self, fn, log_time=True, **args):
        f = io.StringIO()
        with redirect_stdout(f):
            fn(**args)
        self.statement(f.getvalue(), log_time=log_time, echo=False)

    def log(self, phrase):
        """logger something that happened.  The first time phrase is passed the
        start time is saved.  The second time the phrase is logged, the
        elapsed time is written

        Parameters
        ----------
            phrase : str
                the thing that happened
        """
        pass
        t = datetime.now()
        if phrase in self.items.keys():
            s = str(t) + ' finished: ' + str(phrase) + " took: " + \
                str(t - self.items[phrase]) + '\n'
            if self.echo:
                print(s,end='')
            if self.filename:
                self.f.write(s)
                self.f.flush()
            self.items.pop(phrase)
        else:
            s = str(t) + ' starting: ' + str(phrase) + '\n'
            if self.echo:
                print(s,end='')
            if self.filename:
                self.f.write(s)
                self.f.flush()
            self.items[phrase] = copy.deepcopy(t)

    def warn(self, message):
        """write a warning to the logger file.

        Parameters
        ----------
        message : str
            the warning text
        """
        s = str(datetime.now()) + " WARNING: " + message + '\n'
        if self.echo:
            print(s, end='')
        if self.filename:
            self.f.write(s)
            self.f.flush
        warnings.warn(s, RuntimeWarning)

    def lraise(self, message):
        """logger an exception, close the logger file, then raise the exception

        Parameters
        ----------
        message : str
            the exception message

        Raises
        ------
            exception with message
        """
        s = str(datetime.now()) + " ERROR: " + message + '\n'
        print(s,end='')
        if self.filename:
            self.f.write(s)
            self.f.flush
            self.f.close()
        raise Exception(message)
