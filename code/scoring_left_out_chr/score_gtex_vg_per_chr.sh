#!/usr/bin/env bash

hailctl dataproc start qbwcluster -w 40 --max-idle 10m
hailctl dataproc submit qbwcluster ~/PycharmProjects/python3projects/gtex_finemapping/ems_geuvadis_validation/score_gtex_vg_prep_per_chr.py
hailctl dataproc stop qbwcluster

