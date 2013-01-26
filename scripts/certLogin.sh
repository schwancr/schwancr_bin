#!/bin/bash

screen -d -m ssh -A -o ServerAliveInterval=180 -o Ciphers=arcfour -N -L12000:certainty-login.stanford.edu:22 vspm32
