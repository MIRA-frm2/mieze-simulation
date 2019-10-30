import numpy as np
import csv

def transform_cartesian_to_cylindrical(x, y, z):
    """Transform coordinates from cartesian to cylindrical."""
    x = x
    rho = np.sqrt(y ** 2 + z ** 2)
    if z:
        theta = np.arctan(y / z)
    else:
        theta = 0
    return x, rho, theta


def transform_cylindrical_to_cartesian(x, rho, theta):
    """Transform coordinates from cylindrical to cartesian."""
    x = x
    y = rho * np.cos(theta)
    z = rho * np.sin(theta)
    return x, y, z


def save_data_to_file(data, file_name, extension='.csv'):
    full_filename = f'{file_name}{extension}'

    with open(full_filename, 'w') as file:

        print(f'Writing data to file {full_filename}')

        csv_writer = csv.writer(file, delimiter=',')
        csv_writer.writerow(["x", "y", "z", "Bx", "By", "Bz"])

        for point, field in data.items():
            row = list(point) + list(field)
            csv_writer.writerow(row)


def read_data_from_file():
    pass