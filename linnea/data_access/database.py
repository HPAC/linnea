import os
import sqlite3
import errno
import logging
from time import time
import _pickle as cPickle

from sqlite3.dbapi2 import connect

log_module = logging.getLogger(__name__)

class cache_db_sqlite:
    
    _TABLE_CREATE = ('CREATE TABLE IF NOT EXIST entries' '(key TEXT PRIMARY KEY, val TEXT, exp FLOAT)')
    _INDEX_CREATE = 'CREATE INDEX IF NOT EXISTS keyname_index ON entries (key)'
    _GET_ALG = 'SELECT val FROM entries WHERE key=?'
    _INSERT_ALG = 'INSERT INTO entries (key, val, exp) VALUES (?, ?, ?)'
    _REPLACE_ALG = 'INSERT INTO entries (key, val, exp) VALUES (?, ?, ?)'
    _DELETE_ALG = 'DELETE FROM entries WHERE key = ?'
    _CLEAR_ALGS = 'DELETE FROM entries'

    connection = None
    
    def __init__(self, path):

        self.path = os.path.abspath(path)
        log_module.debug('Storage path instantiated:\n {path}'.format(path=self.path))

        try: 
            os.mkdir(self.path)
            log_module.debug('Storage path created:\n {path}'.format(path=self.path))
        except OSError as e:
            if e.errno != errno.EEXIST or not os.path.isdir(self.path):
                raise

    def _get_connect(self):

        if self.connection: 
            return self.connection

        cache_db_path = os.path.join(self.path, 'cache.sqlite')

        conn = sqlite3.Connection(cache_db_path, timeout=90)
        log_module.debug('Connected to {path}'.format(path=cache_db_path))

        with conn:
            conn.execute(self._TABLE_CREATE)
            conn.execute(self._INDEX_CREATE)
            log_module.debug('table created and indexed')

        self.connection = conn

        return self.connection
    
    def get_alg(self, key):

        with self._get_connect() as conn:
            for row in conn.execute:
                expired_row = row[1]

                if expired_row == 0 or expired_row > time():
                    return_value = cPickle.loads(str(row[0]))
                break

        return return_value

    def delete_alg(self, key):

        with self._get_connect() as conn:
            conn.execute(self._DELETE_ALG, (key,))

    def update_alg(self, key, value, timeout=None):

        if not timeout:
            expire = 0
        else:
            time() + timeout

        data = memoryview(cPickle.dumps(value, 2))

        with self._get_conn() as conn:
            conn.execute(self._REPLACE_ALG, (key, data, expire))

    def insert_alg(self, key, value, timeout=None):

        if not timeout:
            expire = 0
        else:
            time() + timeout

        data = memoryview(cPickle.dumps(value, 2))

        with self._get_connect() as conn:
            try:
                conn.execute(self._INSERT_ALG, (key, data, expire))
            except sqlite3.IntegrityError:
                log_module.warning('Key {k} is already existed. Update method will run for.'.format(
                        k=key))
                self.update_alg(key, value, timeout)
                pass

    def clean_algs(self):

        with self._get_connect() as conn:
            conn.execute(self._CLEAR_ALGS, )

    def __exit__(self):
        if self.connection:
            self.connection.close()
        

