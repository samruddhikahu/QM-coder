# MATLAB-code-for-QM-coder
This repository contains MATLAB code for a binary arithmetic QM encoder which can be used to encode 8x8 DCT transformed, quantized image sub-blocks.

Useage of this QM coder is shown in the demo_JPEG_QM.m file.

Input images must be numbered as Image1.bmp, Image2.bmp, ... and so on. 
They must be stored in the folder "D:/MATLAB/JPEG with ArithmeticQM/Input Images".
You may change the location where the input images are stored. 
But make sure that you change the folder in the code accordingly.

Reconstructed output images will be stored in the folder "D:/MATLAB/JPEG with ArithmeticQM/Reconstructed Images".
Quality parameters of the compressed images are stored in an excel file named "Data1.xls" in the same folder.
