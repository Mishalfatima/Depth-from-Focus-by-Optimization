# Depth-from-Focus-by-Optimization
Step 1: Image alignment

1.	The image with the nearest focus is used as first reference. 
2.	Every image is loaded into MATLAB and SURF descriptors are computed using Image Alignment Toolbox. 
3.	We use RANSAC from Image Alignment Toolbox to compute inliers. 
4.	Then, we use inverse warping function from Image Alignment Toolbox to compute warped images. 
5.	All these images are saved separately and will be used in subsequent steps. 

Step 2: Focus Measure
Focus measure means to map an image to a value that represents the degree of focus of the measure. The algorithm was implemented from the paper: “A novel algorithm for estimation of depth map using image focus for 3D shape recovery in the presence of noise”.

1.	First of all, Optical Transfer function is computed. 
2.	Afterwards, the spectrum of each aligned image is found using fft2 command in MATLAB. These images are then element wise multiplied with OTF function found in step 1. 
3.	After element-wise multiplication, inverse Fourier transform is taken by using ifft2 command in MATLAB.
4.	All resulting images are concatenated together to form a 4D vector (row, col, channel, no. of images).
5.	Mp is a 2D matrix representing the maximum value of focus at every pixel location in the concatenated stack of images. Mf is also a 2D matrix representing the frame number in which the maximum focus value is present corresponding to every pixel location.   

Step 3: Graph cuts
 Graph cut optimization can be employed to efficiently solve a wide variety of low-level computer vision problems such as image smoothing, the stereo correspondence problem, image segmentation, and many other computer vision problems that can be formulated in terms of energy minimization.

1.	First step is to find out data and smoothing cost term by making use of Mf 2D matrix as computed in step 2.
2.	Data cost represents the cost incurred from the current frame to the frame with maximum focus value.
3.	Next, we find out smoothing cost.
4.	Once we have both costs, we pass this information to Graph cuts function (in-built library) which returns depth. We use the depth vector to compute an all-in-Focus measure.  

Step 4: All-in-Focus Measure

1.	The result from Graph cuts are used to extract pixel values from the original aligned images corresponding to the frame numbers and form an all-in-focus image. 

Step 5: Depth Refinement

1.	Output from graph cuts is padded and a median filter is passed through the whole image to find refined depth.

References:

1. “A novel algorithm for estimation of depth map using
image focus for 3 D shape recovery in the presence of noise, A Malik
and T S Choi, Pattern Recognition 2008

2.“An Experimental Comparison of Min Cut/Max Flow
Algorithms for Energy Minimization in Vision, Yuri Boykov and Vla
dimir Kolmogorov, IEEE TPAMI 2004”




