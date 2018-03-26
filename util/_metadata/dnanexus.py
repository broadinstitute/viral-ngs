#!/usr/bin/env python

"""Convert dnanexus data to our data format"""

import json
import dxpy.api


def convert_dnanexus_data():
    z = dxpy.api.system_find_executions({
        'class': 'job',
        'project': 'project-F72VZPQ0PPyf3Z54GyF0gkzG',
        'state': 'done',
        'describe': True,
    })['results']

    json_str = json.dumps(z, indent=4)
    print(json_str)

if __name__ == '__main__':
    convert_dnanexus_data()





