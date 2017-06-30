import gzip
import logging
from StringIO import StringIO

import numpy


def detect_schema(filename):
    logging.info('reading gzip file started from ' + filename)
    f = gzip.open(filename, 'rb')
    file_content = f.read()
    f.close()
    logging.info('genfromtxt started')
    narray = numpy.genfromtxt(StringIO(file_content), names=True, delimiter=',', dtype=None)
    logging.info('genfromtxt finished')
    return generate_schema(narray)


def parse_and_convert(filename, global_schema):
    result = []
    first_line = True
    sorted_schema = None
    with gzip.open(filename, 'r') as fin:
        for line in fin:
            parts = line.strip().split(',')
            if first_line:
                sorted_schema = sort_schema(global_schema, parts)
                first_line = False
                continue
            values = {}
            for i in range(0, len(parts)):
                field = sorted_schema[i]
                value = parts[i]
                field_name = field['name']
                field_type = field['type']
                if not value:
                    continue
                values[field_name] = cast_value(value, field_type)
            result.append(values)
    return result


def cast_value(value, field_type):
    if field_type == 'double':
        return float(value)
    if field_type == 'int64':
        return long(value)
    if field_type == 'boolean':
        return bool(value)
    if field_type == 'string':
        return value
    raise Exception('unknown type ' + field_type)


def sort_schema(global_schema, names):
    schema_dict = {}
    for column in global_schema:
        schema_dict[column['name']] = column
    sorted_schema = []
    for name in names:
        sorted_schema.append(schema_dict[name])
    return sorted_schema


def convert(narray):
    array = []
    keys = narray.dtype.names
    logging.info('Converting to array started')
    for row in narray.tolist():
        values = {}
        for i in range(0, len(keys)):
            key = keys[i]
            value = row[i]
            if value == value:
                values[key] = value
            else:
                values[key] = None
        array.append(values)
    logging.info('Converting finished')
    return array


def generate_schema(narray):
    # logging.info('Generating schema')
    schema = []
    for name in narray.dtype.names:
        column = {'name': name}
        t = narray.dtype[name]
        if t == 'float64':
            column['type'] = 'double'
        elif t == 'int':
            column['type'] = 'int64'
        elif t == 'bool':
            column['type'] = 'boolean'
        elif t.kind == 'S':
            column['type'] = 'string'
        else:
            raise Exception('unknown type ' + t)
        schema.append(column)
    return schema


def generate_schema_from_article():
    # https://gaia.esac.esa.int/documentation/GDR1/Catalogue_consolidation/sec_cu1cva/cu9gat_tables.html
    schema = []
    with open('data_fields.txt') as f:
        content = f.readlines()
        for l in content:
            parts = l.strip().split('\t')
            field_name = parts[0]
            field_type = parts[2]
            schema.append({
                'name': field_name,
                'type': convert_field_type_to_yt_format(field_type)
            })
    return schema


def convert_field_type_to_yt_format(field_type):
    if field_type == 'long' or field_type == 'int' or field_type == 'short':
        return 'int64'
    if field_type == 'float':
        return 'double'
    return field_type


if __name__ == '__main__':
    global_schema = generate_schema_from_article()
    result = parse_and_convert('gaia_source/GaiaSource_000-000-000.csv.gz', global_schema)
    print result[0]
    #
    # for i in range(0, 10):
    #     print array[i]
