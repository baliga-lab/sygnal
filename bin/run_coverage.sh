#!/bin/bash

cd python && PYTHONPATH=src coverage3 run test/stats_test.py $@

