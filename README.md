# Differentiated Explanation

This is the pytorch implemention for paper “Differentiated Explanation of Deep Neural Networks with Skewed Distributions”. We propose a simple but efficent approch for the differentiated explanations of black-box classifiers.

To do this, we introduce a trainable relevance estimator that produces relevance scores in a skewed distribution. Specifically, we present the concept of distribution controllers and integrate it with a neural network to directly guide the distribution of relevance scores. By analyzing the effect of the skewness of distributions, we develop two types of controllers with right-skewed distributions for differentiated saliency maps. Then we introduce the classification loss to optimize the estimator. The benefit of this strategy is to better mimic the behavior of deep neural networks without non-trivial hyperparameter tuning, leading to higher faithfulness of explanation.

### Illustration of Differentiated Explanation

<p align="left">
<img src="https://github.com/DEX-DNN/DEX/blob/master/Examples.png" img width="350" height="200" />
</p>

### Benifits of Skewed Distributions
<p align="left">
<img src="https://github.com/DEX-DNN/DEX/blob/master/Distributions.png" img width="380" height="220" />
</p>

### Framework
<p align="left">
<img src="https://github.com/DEX-DNN/DEX/blob/master/Framework.png" img width="400" height="280" />
</p>

### Running the Demo
Our code is implemented on Pytorch 1.1.0. 
```
python saliency_demo.py
```
