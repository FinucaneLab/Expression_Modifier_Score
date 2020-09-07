#!/usr/bin/env bash

hailctl dataproc start qbwcluster -w 40 --max-idle 10m
hailctl dataproc submit qbwcluster ~/PycharmProjects/python3projects/gtex_finemapping/ems_pipe_test_20200125/score_gtex_vg_prep.py
