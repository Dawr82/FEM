from math import sqrt


# For purposes of Gaussian integration.
SCHEMA = {
    1 :  [[1, 1],[-1/sqrt(3), 1/sqrt(3)]],
    2 :  [[5/9, 8/9, 5/9],[-sqrt(3/5), 0, sqrt(3/5)]],
    3 :  [[0.347855, 0.652145, 0.347855, 0.652145],[-0.861136, -0.339981, 0.339981, 0.0861136]],
}


# Gaussian integration - 1D, normalized (x from -1 to 1) coordinate system.
def gaussian_integration1D(f, N):
    result = 0

    try:
        for i in range(N + 1):
            result += f(SCHEMA[N][1][i]) * SCHEMA[N][0][i]
    except KeyError:
        print("No data in SCHEMA for that N.")
   
    return result


# Gaussian integration - 2D, normalized (x from -1 to 1) coordinate system.
def gaussian_integration2D(f, N):
    result = 0

    try:
        for i in range(N + 1):
            for j in range(N + 1):
                result += f(SCHEMA[N][1][i], SCHEMA[N][1][j]) * SCHEMA[N][0][i] * SCHEMA[N][0][j]
    except KeyError:
        print("No data in SCHEMA for that N.")

    return result


