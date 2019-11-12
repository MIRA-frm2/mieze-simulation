# -*- coding: utf-8 -*-
#
# This file is part of MIEZE simulation.
# Copyright (C) 2019 TUM FRM2 E21 Research Group.
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.

"""Plotting scripts from data file."""

import csv


def read_data(file_name='../data/data.csv'):
    x = list()
    y = list()
    z = list()

    bx = list()
    by = list()
    bz = list()

    with open(file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                x.append(row[0])
                y.append(row[1])
                z.append(row[2])

                bx.append(row[3])
                by.append(row[4])
                bz.append(row[5])

    return x, y, z, bx, by, bz


if __name__ == "__main__":
    print(read_data())
