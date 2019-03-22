#!/bin/bash

awk 'function human(x) {
         s=" B   KiB MiB GiB TiB EiB PiB ZiB YiB"
         while (x>=1024 && length(s)>1) 
               {x/=1024; s=substr(s,5)}
         s=substr(s,1,4)
         xf=(s==" B  ")?"%5d   ":"%8.2f"
         return sprintf( xf"%s\n", x, s)
      }
      {gsub(/[0-9]+/, human($1)); print}'
