#!/bin/bash

screen -S vspcert -d -m ssh -o ServerAliveInterval=180 -o Ciphers=arcfour -N -L11000:localhost:12000 schwancr@vspm24
