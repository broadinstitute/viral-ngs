#!/bin/bash

head -n1 /proc/meminfo | awk '{print int($2*0.5/1024)}'
