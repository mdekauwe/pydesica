# pydesica: a simple simulation model for the desiccation dynamics of plants coded in python

The desica model is a simple implementation of a soil-plant-atmosphere water transport model, specifically designed to model the dynamics of the desiccation of plants to mortality.

When soil water uptake does not meet demand by the plant (which includes an important minimum conductance term), water is drawn from reserves held in the stem and leaves. When the water potential of the stem is low enough, the plant dies.

## Key References:

1. Tuzet et al. (2003) A coupled model of stomatal conductance, photosynthesis and transpiration. Plant, Cell and Environment 26, 1097–1116.
2. Xu X, Medvigy D, Powers JS, Becknell JM, Guan K (2016) Diversity in plant hydraulic traits explains seasonal and inter-annual variations of vegetation dynamics in seasonally dry tropical forests. New Phytologist, 212, 80–95.

## Running the model

```bash
$ python src/desica.py
```
