# [posest](https://users.ics.forth.gr/~lourakis/posest/) - *Robust 6DoF Pose Estimation from 3D-2D Correspondences*

This is *posest*, a library for 3D pose estimation from point correspondences that is distributed under the GPLv2.
Pose estimation refers to the computation of position and orientation estimates that fully define the posture
of a rigid object in space (6 DoF in total).
The estimation is based on a set of known 3D points and their corresponding 2D projections on an image.

The library estimates the relative motion between the 3D points and the camera system. Single or binocular camera
systems are supported. Image points typically originate from local features, however posest is oblivious to their origin.
Mismatched 3D-2D point pairs (i.e. outliers) are handled with [RANSAC](https://en.wikipedia.org/wiki/Random_sample_consensus).
More details can be found at [posest's web site](https://users.ics.forth.gr/~lourakis/posest/).

## Required libraries
posest requires the [levmar](https://users.ics.forth.gr/~lourakis/levmar/) library to build.

## Cite as
If you use this code in your published work, please cite the following paper:<br><br>
<pre>
@incollection{Lourakis13model,
  author={Lourakis, Manolis and Zabulis, Xenophon}
  title="{Model-Based Pose Estimation for Rigid Objects}",
  booktitle={Computer Vision Systems},
  isbn={978-3-642-39401-0},
  volume={7963},
  series={Lecture Notes in Computer Science},
  editor={Chen, Mei and Leibe, Bastian and Neumann, Bernd},
  publisher={Springer Berlin Heidelberg},
  doi={10.1007/978-3-642-39402-7_9},
  pages={83--92},
  year="2013",
  url="http://users.ics.forth.gr/~lourakis/posest/ICVS13.pdf"
}
<pre>
