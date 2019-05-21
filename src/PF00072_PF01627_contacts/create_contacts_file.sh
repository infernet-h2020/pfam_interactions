#!/bin/bash

awk '$7<8 {print FILENAME, $0}' * > sc_sc_contacts.txt
