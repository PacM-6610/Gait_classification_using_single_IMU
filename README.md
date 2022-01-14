# Gait_classification_using_single_IMU
This paper explores two classification approaches to detect subtle gaitchanges while walking on a flat plane and when walking on stairs. This is doneby using machine learning and an inertial sensor worn close to the body’s CoM(lower back/ waist). Walking on stairs represents a challenging bio-mechanicaleffort compared to walking on a flat surface. Its detection can help preventfall injuries in people with physical impairment by identifying their ability towalk on stairs and potential risk areas. The first approach takes advantageof feature selection using three traditional filter methods; T-test (TT), Maxi-mum Relevance Minimum Redundancy (MRMR) and chi square test (CST).The second approach takes advantage of feature reduction using a dimension-ality reduction algorithm; principal component analysis (PCA). The resultingmodels of the two approaches were trained and tested using a single inertialmeasurement unit worn in the lower back from the dataset of Luo and al. and the K-folds and leave-one-out cross-validation algorithms