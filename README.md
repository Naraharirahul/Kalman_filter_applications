# Kalman_filter_applications

This repository contains UKF and EKF implementations.

## Unscented Kalman Filter

The unscented Kalman filter (UKF) uses a deterministic sampling technique known as the unscented transformation (UT) to pick a minimal set of sample points (called sigma points) around the mean. The sigma points are then propagated through the non‐linear functions, from which a new mean and covariance estimate are then formed.

### Project

This project implements an Unscented Kalman Filter (UKF) to estimate the state of multiple cars on a highway using noisy lidar and radar measurements.

## Extended Kalman Filter

In the EKF, the state transition and the observation models don't need to be linear functions of the state but may instead be differential functions. Extended Kalman Filter can handle non‐linear relationships; but we still need Gaussian noise model.

### Project

Applied Extended Kalman filter algorithm for tracking a bicycle's position and velocity in a simulator by fusing Lidar sensor
measurements for the update step and radar sensor for measurement step and linearizing the EKF algorithm's data.
