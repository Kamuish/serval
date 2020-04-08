from .read_funcs import read_harps_ccf, read_template


def read_handler(read_type, filename):
    map_inputs = {
        'harps_ccf': read_harps_ccf,
        'template': read_template
    }

    return map_inputs[read_type](filename)