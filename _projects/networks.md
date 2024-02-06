---
layout: page
title: network analysis
description: Developing and applying tools for network analysis in high-dimensional -omics data.  
img: assets/img/circleplot.png
importance: 2
category: work
related_publications: true
---

Network analysis in genomic data is one of my primary research interests. Biological systems are highly complex and intrinsically interconnected, and there is much we can learn about diseases and other conditions by interrogating genomic networks. Many of my research projects have involved applying existing network analysis tools, but I have also worked on developing such tools like [hdWGCNA](https://smorabit.github.io/hdWGCNA/).

***

### Method development

<div class="row justify-content-sm-center">
    <div class="col-sm-8 mt-3 mt-md-0">
One of my primary research interests is in method development for network analysis in genomics data. I am the lead developer of <a href="https://github.com/smorabit/hdWGCNA">hdWGCNA</a>, an R package for co-expression network analysis in single-cell and spatial transcriptomics data. This method constructs networks in specific contexts that reflect the underlying biology of individual cell types or disease conditions. hdWGCNA also includes downstream tools for data visualization and statistical testing. We developed this package to be user-friendly and we leveraged common data formats which enable a borad range of researchers to use this package. To date hdWGCNA has already gathered hundreds of users across many different fields of research. Check out our publication to learn more about this method {% cite Morabito2023 %}.
    </div>
    <div class="col-sm-4 mt-3 mt-md-0">
        {% include figure.liquid path="assets/img/logo_v4.png" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
</div>

***

### Applications

<div class="row justify-content-sm-center">
  <div class="col-sm-4 mt-3 mt-md-0">
        {% include figure.liquid path="assets/img/circleplot.png" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
    <div class="col-sm-8 mt-3 mt-md-0">  
Nearly all of my research projects have involved applying network analysis approaches in some capacity. These approaches enable us to identify systems-level changes in disease or in specific cell states, beyond what differential expression testing is capable of revealing. For example, in our recent preprint {% cite Miyoshi2023 %} we performed multi-scale network analysis across different cortical layers in our spatial dataset to reveal an altered neuroinflammatory signature in disease. In some cases, I have also worked on other types of network analysis aside from co-expression networks, like transcription factor regulatory networks in {% cite Morabito2021 %}. 
    </div>
</div>

***
