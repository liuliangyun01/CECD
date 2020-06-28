CECD

The soft code called CECD (Classification Extension-based Cloud Detection method ) is used for automated clouds masking for Landsat-8, Sentinel-2 or other medium-resolution satellite images (Such as GaoFens).

CECD is to extend the local adaptive classifier developed at a low-resolution level to the corresponding medium-resolution image for cloud masking. 

If you have any questions, please contact Liangyun Liu (liuly@radi.ac.cn) and Xidong Chen (chenxd@radi.ac.cn) at Aerospace Information Research Institute, Chinese Academy of Sciences.

IMPORTANT:

This Github page ONLY includes the IDL code for CECD V1.0. Note that, because IDL cannot save and export the trained model, the code here only includes the process of training sample generation. The final classification process needs to be conducted manually based on the selected samples. 

The input data required by our program are: 300mPRoba-V TOA reflectance data and resampled 300m Landsat-8/Sentinel-2 (TOA) reflectance data. The output data are THOT image, THMFCM image (THOT image after MFCM correction), and a .roi file that stores the collected training samples.
Parameter:poss1 and poss2 represents the band index of blue, red and SWIR bands in Proba-V and Landsat-8/sentinel-2 images, respectively.

A pair of test images (test_data.rar) are provided in this Github page for code debugging. Two Tif images of THOT and THMFCM are provided to illustrate the effect of THOT and MFCM.  And a Random Forest model (l8_model.pkl) trained by using our collection samples and Landsat-8 images is also provided, and it can be loaded in Python as: with open("l8_model.pkl", "rb") as f: model = pickle.load(f). The input features for the model are: coastal, blue, green, red, nir, swir1, swir2, cirrus. 

Due to the urgent time, this code is only an early version and not user-friendly. In the future, we will update this version of the code. 
