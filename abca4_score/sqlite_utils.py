"""SQLite counterpart of the original project's utils/mysql.py.

Only the functions actually used by the scoring pipeline are ported:
connect / search_db / error_intolerant_search / hard_landing_search /
get_column_names / store_or_update. The return conventions of the MySQL
originals are preserved (search_db and error_intolerant_search return False
when a query has no result rows; hard_landing_search raises instead).
"""

import sqlite3


def connect(db_path: str):
    """Open the SQLite database and return (db, cursor), autocommitting like
    the original abca4_connect() did (it ran `set autocommit=1`)."""
    db = sqlite3.connect(db_path, autocommit=True)
    cursor = db.cursor()
    return db, cursor


def search_db(cursor, qry, verbose=False):
    try:
        cursor.execute(qry)
    except sqlite3.Error as e:
        raise Exception(f"Error running cursor.execute() for qry:\n{qry}\n{e}")

    rows = cursor.fetchall()
    if len(rows) == 0:
        if verbose:
            print(f"No return for query:\n{qry}")
        return False
    return [list(row) for row in rows]


def error_intolerant_search(cursor, qry):
    # in the sqlite3 module errors surface as exceptions already,
    # so this is just search_db with the original name kept for the callers
    return search_db(cursor, qry)


def hard_landing_search(cursor, qry):
    ret = search_db(cursor, qry)
    if not ret:
        raise Exception(f"Hard Landing Search failed (no rows):\n{qry}")
    return ret


def get_column_names(cursor, table_name):
    rows = hard_landing_search(cursor, f"pragma table_info({quote_identifier(table_name)})")
    return [row[1] for row in rows]


def quote_identifier(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def store_or_update(cursor, table, fixed_fields, update_fields, verbose=False,
                    primary_key='id', fail_if_no_exist=False) -> list | bool:
    """Insert a row described by fixed_fields + update_fields, or, if a row
    matching fixed_fields exists, update it with update_fields.
    Returns the primary key(s) of the affected row(s). Same interface as the
    MySQL original, but implemented with parametrized queries."""
    conditions = []
    values = []
    for k, v in fixed_fields.items():
        if v is None:
            conditions.append(f"{quote_identifier(k)} is null")
        else:
            conditions.append(f"{quote_identifier(k)} = ?")
            values.append(v)

    qry = f"select {quote_identifier(primary_key)} from {table}"
    if conditions:
        qry += " where " + " and ".join(conditions)
    if verbose: print(qry, values)
    cursor.execute(qry, values)
    rows = cursor.fetchall()

    if rows:  # exists: update
        primary_keys = [row[0] for row in rows]
        if not update_fields:
            if verbose: print("exists; no update requested")
            return primary_keys
        if verbose: print("exists; updating")
        assignments = ", ".join(f"{quote_identifier(k)} = ?" for k in update_fields)
        placeholders = ", ".join("?" for _ in primary_keys)
        qry = f"update {table} set {assignments} where {quote_identifier(primary_key)} in ({placeholders})"
        cursor.execute(qry, list(update_fields.values()) + primary_keys)
        return primary_keys

    # does not exist: insert
    if fail_if_no_exist:
        raise Exception(f"{fixed_fields} not found in the database")
    if verbose: print("does not exist; making new one")
    fields = dict(fixed_fields)
    fields.update(update_fields or {})
    columns = ", ".join(quote_identifier(k) for k in fields)
    placeholders = ", ".join("?" for _ in fields)
    qry = f"insert into {table} ({columns}) values ({placeholders})"
    cursor.execute(qry, list(fields.values()))
    return [cursor.lastrowid]
