# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [UNRELEASED] - 2024-??-??

## [0.1.1] - 2024-02-29

* Fixed bug in `ComputeLayout` when normalizing or projecting layout coordinates
* Added layout utility functions:
  * `center_layout_coordinates` - Center layout coordinates at the origin
  * `normalize_layout_coordinates` - Normalize layout coordinates such that they are centered at origin and the have a median length (euclidean norm) of 1
  * `project_layout_coordinates_on_unit_sphere` - Project layout 3D coordinates onto the unit sphere

## [0.1.0] - 2024-02-21

* First public release of pixelatorR.
