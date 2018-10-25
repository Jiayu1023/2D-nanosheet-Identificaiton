##########Denoise Raw Image##########
import cv2 as cv
import numpy as np
from matplotlib import pyplot as plt

src=cv.imread("F:/Calculation Results/1.jpg")
cv.namedWindow("input image",cv.WINDOW_AUTOSIZE)
cv.imshow("input image",src)

Res=cv.bilateralFilter(src,0,20,20)
cv.imshow("bi_demo",Res)
cv.waitKey(0)
cv.imwrite("F:/Calculation Results/1 EPF.jpg",Res)
