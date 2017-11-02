import os
from lsst.sims.catUtils.utils import AlertProcessor

processor = AlertProcessor()

schema_name = os.path.join("/Users", "danielsf", "physics", "lsst_171025",
                           "Development", "sample-avro-alert",
                           "schema", "diasource.avsc")

if not os.path.exists(schema_name):
    raise RuntimeError("%s does not exist" % schema_name)

processor.load_diasource_schema(schema_name)

processor.process('test_data_0', 'avro_out_dir/first')
