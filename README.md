# Efficient motion planning
To guarantee safe and efficient driving for automated vehicles in complicated traffic conditions, the motion planning module of automated vehicles are expected to generate collision-free driving policies as soon as possible in varying traffic environment. However, there always exist a tradeoff between efficiency and accuracy for the motion planning algorithms. Besides, most motion planning methods cannot find the desired trajectory under extreme scenarios (e.g., lane change in crowded traffic scenarios). This study proposed an efficient motion planning strategy for automated lane change based on Mixed-Integer Quadratic Optimization (MIQP) and Neural Networks. We modeled the lane change task as a mixed-integer quadratic optimization problem with logical constraints, which allows the planning module to generate feasible, safe and comfortable driving actions for lane changing process. Then, a hierarchical machine learning structure that consists of SVM-based classification layer and NN-based action learning layer is established to generate desired driving policies that can make online, fast and generalized motion planning. Our model is validated in crowded lane change scenarios through numerical simulations and results indicate that our model can provide optimal and efficient motion planning for automated vehicles


In ACC file, we didn’t consider the interaction between each vehicle because as in most lane change decision algorithm, they ignore the interaction to simplify the computation process. And for the surrounding vehicles, we assume them to drive in a constant acceleration. 

In IDM file, we consider the interaction between vehicles, i.e. under this assumption, the obstacle vehicles are no longer ‘blind’ which can take reaction to the dynamic environment. This indicates that our method is worth of further exploration by using more realistic and interactive traffic flow model.


# Refered Paper

\url{https://arxiv.org/abs/1904.08784}
