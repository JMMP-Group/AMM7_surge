.. todo::



NEMO_ for *Nucleus for European Modelling of the Ocean* is a state-of-the-art modelling framework for
research activities and forecasting services in ocean and climate sciences,
developed in a sustainable way by a European consortium since 2008.

.. contents::
   :local:

Overview
========

The NEMO ocean model has 3 major components:

- |OCE| models the ocean {thermo}dynamics and solves the primitive equations
  (:file:`./src/OCE`)
- |ICE| simulates sea-ice {thermo}dynamics, brine inclusions and
  subgrid-scale thickness variations (:file:`./src/ICE`)
- |MBG| models the {on,off}line oceanic tracers transport and biogeochemical processes
  (:file:`./src/TOP`)

These physical core engines are described in
their respective `reference publications <#project-documentation>`_ that
must be cited for any work related to their use (see :doc:`cite`).

Assets and solutions
====================

Not only does the NEMO framework model the ocean circulation,
it offers various features to enable

- Create :doc:`embedded zooms<zooms>` seamlessly thanks to 2-way nesting package AGRIF_.
- Opportunity to integrate an :doc:`external biogeochemistry model<tracers>`
- Versatile :doc:`data assimilation<da>`
- Generation of :doc:`diagnostics<diags>` through effective XIOS_ system
- Roll-out Earth system modeling with :doc:`coupling interface<cplg>` based on OASIS_

Several :doc:`built-in configurations<cfgs>` are provided to
evaluate the skills and performances of the model which
can be used as templates for setting up a new configurations (:file:`./cfgs`).

The user can also checkout available :doc:`idealized test cases<tests>` that
address specific physical processes (:file:`./tests`).

A set of :doc:`utilities <tools>` is also provided to {pre,post}process your data (:file:`./tools`).

Project documentation
=====================

A walkthrough tutorial illustrates how to get code dependencies, compile and execute NEMO
(:file:`./INSTALL.rst`).

Reference manuals and quick start guide can be build from source and
exported to HTML or PDF formats (:file:`./doc`) or
downloaded directly from the :forge:`development platform<wiki/Documentations>`.

============ ================== ===================
 Component    Reference Manual   Quick Start Guide
============ ================== ===================
 |NEMO-OCE|   |DOI man OCE|_     |DOI qsg|
 |NEMO-ICE|   |DOI man ICE|
 |NEMO-MBG|   |DOI man MBG|
============ ================== ===================

Since 2014 the project has a `Special Issue`_ in the open-access journal
Geoscientific Model Development (GMD) from the European Geosciences Union (EGU_).
The main scope is to collect relevant manuscripts covering various topics and
to provide a single portal to assess the model potential and evolution.

Used by a wide audience,
numerous :website:`associated projects<projects>` have been carried out and
extensive :website:`bibliography<bibliography/publications>` published.

Development board
=================

The NEMO Consortium pulling together 5 European institutes
(CMCC_, CNRS_, MOI_, `Met Office`_ and NERC_) plans the sustainable development in order to
keep a reliable evolving framework since 2008.

It defines the |DOI dev stgy|_ that is implemented by the System Team on a yearly basis
in order to release a new version almost every four years.

When the need arises, :forge:`working groups<wiki/WorkingGroups>` are created or resumed to
gather the community expertise for advising on the development activities.

.. |DOI dev stgy| replace:: multi-year development strategy

Disclaimer
==========

The NEMO source code is freely available and distributed under
:download:`CeCILL v2.0 license <../../../LICENSE>` (GNU GPL compatible).

You can use, modify and/or redistribute the software under its terms,
but users are provided only with a limited warranty and the software's authors and
the successive licensor's have only limited liability.
