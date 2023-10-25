# tnfm

### Input File

column | input variable       | unit
-------|----------------------|-----
1      | world time           | s
2      | model time           | s
3      | surface temperature  | K
4      | surface density      | kg/m3
5      | solid accumulation   | m weq./s
6      | liquid accumulation  | m weq./s
7      | surface grain radius | m


### Init File

column | init variable        | unit
-------|----------------------|------
1      | depth                | m
2      | density              | kg/m3
3      | temperature          | K
4      | grain radius         | m
5      | liquid water content | 1
6      | age                  | s

The first line is the base layer, the last line is the top layer of the firn column.


### Output File

column | output variable     | unit
-------|---------------------|-----
1      | depth               | m
2      | density             | kg/m3
3      | temperature         | K
4      | grain radius        | m
5      | age                 | s
6      | luqid water content | 1
