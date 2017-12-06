#!/bin/bash

head -n1 /proc/meminfo | awk '{print int($2*0.9/1024/1024)}'
