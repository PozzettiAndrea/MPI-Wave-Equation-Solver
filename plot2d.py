import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import imageio
import sys

num = int(sys.argv[1])
t_out = 0.25

row = int(os.listdir("./out")[0][5])
col = int(os.listdir("./out")[0][7])

for i in range(1,num):
    grid = pd.DataFrame()
    for k in range(row):
        line = pd.DataFrame()
        for j in range(col):
            square = pd.read_csv('./out/id_%d_%dx%d_output_%d.dat' % ((k*row + j), row, col, i), header=None, sep="\t", dtype=float)
            square = square.drop([np.shape(square)[1]-1], axis=1)
            line = pd.concat((line, square), axis=1, ignore_index=True)
        grid = pd.concat((grid, line), ignore_index=True)
    
    #ADAPTED FROM STACKOVERFLOW TOO
    fig, ax = plt.subplots()
    ax.imshow(grid ,cmap='Blues')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('T: %.1f' % (t_out*i),fontsize='large')
    plt.savefig("./images/foo" + ("0" * (4-len(str(i)))) + str(i) + ".png")
    plt.close(fig)

#I HAVE COPIED THIS FROM STACKOVERFLOW https://stackoverflow.com/questions/41228209/making-gif-from-images-using-imageio-in-python
png_dir = 'images/'
images = []
for file_name in sorted(os.listdir(png_dir)):
    if file_name.endswith('.png'):
        file_path = os.path.join(png_dir, file_name)
        images.append(imageio.imread(file_path))
imageio.mimsave('images/animation.gif', images)