# WCGAN-GP-Imagenes-Coloidales
The Wasserstein Conditional Generative Adversarial Network with Gradient Penalty, abbreviated as WCGAN-GP, is a model developed by  [Walia, Tierney and McKeever 2020](https://www.researchgate.net/publication/347437993_Synthesising_Tabular_Data_using_Wasserstein_Conditional_GANs_with_Gradient_Penalty_WCGAN-GP)   for data generation. The WCGAN-GP uses the Wasserstein loss along with the gradient penalty. This approach contributes to greater stability during training. In addition, it stands out for being a conditional GAN, which implies its ability to generate conditional data according to the input label.

This repository hosts the implementation of WCGAN-GP and exemplifies its application by generating synthetic data from a specific data set (Images obtained through bright field video microscopy). It is worth noting that, in this context, WCGAN-GP has been adapted to generate statistically equivalent quasi-2D colloidal images.

