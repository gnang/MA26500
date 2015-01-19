︠965f2f88-c9b0-4b96-9027-c9b8ae64ec33ai︠
%hide
%auto
DATA="foo.data/"
︠dd956472-9200-42e7-b6b3-1e12fc7ce8fa︠
# Importing image pocessing libraries
import pylab, numpy
from scipy import misc

ImPi = pylab.imread('foo.data/pi_icon.png')
ImHd = pylab.imread('foo.data/head_icon.png')
g = graphics_array([matrix_plot(ImPi),matrix_plot(ImHd)])
g.show(figsize=[5,5])
︠ed086b32-833d-44d8-a7b2-9d0ddacd120c︠

︠240d9825-954f-411a-9b41-a71815b533f1︠
matrix_plot(0.5*ImPi+0.5*ImHd)
︠51681ac4-1d60-4687-8993-d9535d4dedc9︠
︠98cae98a-aad0-4d00-840f-f4b7b30fc2f7︠
X = pylab.imread('foo.data/smile.png')
matrix_plot(X)
︠5ac3f91b-f154-4f1b-b6dc-f6c8385b668b︠
# Initializing the matrix
Id15 = identity_matrix(15)
Rs = Id15.tensor_product(hadamard_matrix(2))

# Permuting the column of the matrix
for i in range(1,15):
    tmp = Rs[:,i]
    Rs[:,i] = Rs[:,2*i]
    Rs[:,2*i] = tmp
    
# Ask not what the product of matrices can do for you
# but ask for what you can do with the matrix product !!!
# Ploting the result
Y = numpy.zeros(X.shape)
for i in range(3):
    Y[:,:,i] = 0.25*numpy.dot(numpy.dot(Rs.transpose(),X[:,:,i]),Rs)

# Ploting the resulting matrix
matrix_plot(Y)









