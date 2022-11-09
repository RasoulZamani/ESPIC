"""
it is for testing things in dev proccess.
it will be deleted after deploy!
"""
import os
from parameters import *
from sys import platform
import time
from multiprocessing.connection import wait
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
import os

x =[5, 7, 8, 7, 2, 17, 2, 9,
	4, 11, 12, 9, 6]

y =[99, 86, 87, 88, 100, 86,
	103, 87, 94, 78, 77, 85, 86]


"""
for n in range(4):
    fig, ax = plt.subplots()
    ax.scatter(x,y)
    ax.set( title=f" fig {n}" )
    plt.show()
    time.sleep(1)

def plt_xv(n):
    fig = plt.figure()
    plt.xlabel('Distance')
    plt.ylabel('Elevation')
    plt.scatter(x,y)

    #fig.show()
    fig.savefig(r"results/fig"+f"_{n}")

    
dir="./results"
os.makedirs(dir, exist_ok=True)  
files = os.listdir(dir)
if len(files)>0:
    for f in files:
        file_adrs = os.path.join(dir,f)
        os.remove(file_adrs)

a= os.path.join(RESULT_DIR)
print(a)

""" 
"""
class Person():
    def __init__(self,name,age):
        self.name = name
        self.age= age
    def display(self):
        print(self.name, self.age)

person = Person(name="rasoul", age=30)
person.display()

class Student(Person):
    def __init__(self,name,age,score):
        super().__init__(name, age)
        self.score=score
    def display(self):
        print(self.score)

stud = Student(name="ali", age= 29, score=20)
stud.display()
"""  

print(2.1 ** (0.5))