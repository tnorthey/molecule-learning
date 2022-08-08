import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression

# read Coulomb matrices as x
# ...
# read IAM vectors as y
# ...

x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=4, random_state=4)


