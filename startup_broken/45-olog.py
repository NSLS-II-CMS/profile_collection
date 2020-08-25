


simple_template = """{{- start.plan_name }} ['{{ start.uid[:6] }}'] (scan num: {{ start.scan_id }})"""


count_template = """{% if 'plan_header_override' in start -%}{{start.plan_header_override}}{% else %}{{start.plan_name}}{% endif %} {% if 'sample_savename' in start -%}{{start.sample_savename}}{% else %}{%if 'sample_name' in start -%}{{start.sample_name}}{% endif %}{% endif %} :  {{start.plan_args.num}} ['{{ start.uid[:6] }}'] (scan num: {{ start.scan_id }})

Scan Plan
----------------------------------------
{{ start.plan_type }}
{%- for k, v in start.plan_args | dictsort %}
    {{ k }}: {{ v }}
{%-  endfor %}
{% if 'signature' in start -%}
Call:
    {{ start.signature }}
{% endif %}

Metadata
----------------------------------------
{% for k, v in start.items() -%}
{%- if k not in ['plan_type', 'plan_args'] -%}{{ k }} : {{ v }}
{% endif -%}
{%- endfor %}



"""



#single_motor_template = """{{- start.plan_name}} :  {{ start.motors[0]}} {{start.plan_args.start}} {{start.plan_args.stop}} {{start.plan_args.num}} ['{{ start.uid[:6] }}'] (scan num: {{ start.scan_id }})
single_motor_template = """{% if 'plan_header_override' in start -%}{{start.plan_header_override}}{% else %}{{start.plan_name}}{% endif %} :  {{ start.motors[0]}}  {{'%0.3f' %start.plan_args.start|float}}    {{'%0.3f' %start.plan_args.stop|float}} {{start.plan_args.num}} ['{{ start.uid[:6] }}'] (scan num: {{ start.scan_id }})


Scan Plan
----------------------------------------
{{ start.plan_type }}
{%- for k, v in start.plan_args | dictsort %}
    {{ k }}: {{ v }}
{%-  endfor %}
{% if 'signature' in start -%}
Call:
    {{ start.signature }}
{% endif %}

Metadata
----------------------------------------
{% for k, v in start.items() -%}
{%- if k not in ['plan_type', 'plan_args'] -%}{{ k }} : {{ v }}
{% endif -%}


{%- endfor -%}

"""


from collections import defaultdict
TEMPLATES = defaultdict(lambda: simple_template)
TEMPLATES['ct'] = count_template
TEMPLATES['count'] = count_template
TEMPLATES['scan'] = single_motor_template
TEMPLATES['dscan'] = single_motor_template
TEMPLATES['ascan'] = single_motor_template
TEMPLATES['ID_calibration'] = single_motor_template

#from jinja2 import Template # What is this needed for?

# connect olog
from functools import partial
from pyOlog import SimpleOlogClient
from bluesky.callbacks.olog import logbook_cb_factory

# Set up the logbook. This configures bluesky's summaries of
# data acquisition (scan type, ID, etc.).

LOGBOOKS = ['Data Acquisition']  # list of logbook names to publish to
simple_olog_client = SimpleOlogClient()
generic_logbook_func = simple_olog_client.log
configured_logbook_func = partial(generic_logbook_func, logbooks=LOGBOOKS)

# This is for ophyd.commands.get_logbook, which simply looks for
# a variable called 'logbook' in the global IPython namespace.
logbook = simple_olog_client


#logbook_cb = logbook_cb_factory(configured_logbook_func)
logbook_cb = logbook_cb_factory(configured_logbook_func, desc_dispatch=TEMPLATES)


# Note: According to Yugang, to get the filename from the databroker (after the end of a run), you can do something like:
#header = db[-1] # Get most recent record
#events = db.get_events(header)
#  events_list = [ ev for ev in events ]
#keys = [k for k, v in header.descriptors[0]['data_keys'].items()     if 'external' in v]
#key, = keys
#filenames =  [  str( ev['data'][key][0]) + '_'+ str(ev['data'][key][2]['seq_id']) for ev in events]     
# Refer to Yugang's chxanalys package for full code.
#    https://github.com/yugangzhang/chxanalys
from functools import partial
from pyOlog import SimpleOlogClient
import queue
import threading
from warnings import warn


def submit_to_olog(queue, cb):
    while True:
        name, doc = queue.get()  # waits until document is available
        try:
            cb(name, doc)
        except Exception as exc:
            warn('This olog is giving errors. This will not be logged.'
                 'Error:' + str(exc))

olog_queue = queue.Queue(maxsize=100)
olog_thread = threading.Thread(target=submit_to_olog, args=(olog_queue, logbook_cb), daemon=True)
olog_thread.start()

def send_to_olog_queue(name, doc):
    try:
        olog_queue.put((name, doc), block=False)
    except queue.Full:
        warn('The olog queue is full. This will not be logged.')

RE.subscribe(send_to_olog_queue, 'start')
