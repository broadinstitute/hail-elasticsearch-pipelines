#!/usr/bin/env bash

set -x

# must run sudo /sbin/sysctl -w vm.max_map_count=262144  on the VM

mkdir -p /logs
chown elasticsearch /elasticsearch-data /logs

su elasticsearch -c "/usr/local/elasticsearch-${ELASTICSEARCH_VERSION}/bin/elasticsearch \
    -E network.host=0.0.0.0 \
    -E http.port=${ELASTICSEARCH_SERVICE_PORT} \
    -E path.data=/elasticsearch-data \
    -E path.logs=/logs"

echo elasticsearch started on port ${ELASTICSEARCH_PORT}!!