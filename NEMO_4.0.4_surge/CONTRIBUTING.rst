************
Contributing
************

.. todo::



.. contents::
   :local:

Sending feedbacks
=================

|  Sending feedbacks is a useful way to contribute to NEMO efficency and reliability. Before doing so,
   please check here :forge:`search <search on the developement platform>` in wiki, tickets, forum and online
   documentation if the subject has already been discussed. You can either contribute to an existing
   discussion, or
| Create an entry for the discussion online, according to your needs

- You have a question: create a topic in the appropriate :forge:`discussion <forum>`
- You would like to raise and issue: open a new ticket of the right type depending of its severity

  - "Unavoidable" :forge:`newticket?type=Bug       <bug>`

  - "Workable"    :forge:`newticket?type=Defect <defect>`

Please follow the guidelines and try to be as specific as possible in the ticket description.

New development
===============

You have build a development relevant for NEMO shared reference: an addition of the source code,
a full fork of the reference, ...

You may want to share it with the community (see Hack below) or to propose it for implementation in the future
NEMO release (see Proposal / Task below).

The proposals for developments to be included in the shared NEMO reference are first examined by NEMO Developers
Committee / Scientific Advisory Board.
The implementation of a new development requires some additionnal work from the intial developer.
These tasks will need to be scheduled with NEMO System Team.


Hack
----

You only would like to inform NEMO community about your developments.
You can promote your work on NEMO forum gathering  the contributions fromof the community by creating
a specific topic here :forge:`discussion/forum/5 <dedicated forum>`


Proposal / Task
---------------

| Your development is quite small, and you would only like to offer it as a possible enhancement. Please suggest it
  as an enhancement here :forge:`newticket?type=Enhancement <enhancement>` . It will be taken in account, if
  feasible, by NEMO System Team. To ease the process, it is suggested, rather than attaching the modified
  routines to the ticket, to highlight the proposed changes by adding to the ticket the output of ``svn diff``
  or ``svn patch`` from your working copy.

| Your development seems relevant for addition into the future release of NEMO shared reference.
  Implementing it into NEMO shared reference following the usual quality control will require some additionnal work
  from you and also from the NEMO System Team in charge of NEMO development. In order to evaluate the work,
  your suggestion should be send as a proposed enhancement here :forge:`newticket?type=Enhancement <enhancement>`
  including description of the development, its implementation, and the existing validations.

  The proposed enhancement will be examined  by NEMO Developers Committee / Scientific Advisory Board.
  Once approved by the  Committee, the assicated development task can be scheduled in NEMO development work plan,
  and tasks distributed between you as initial developer and PI of this development action, and the NEMO System Team.

  Once sucessful (meeting the usual quality control steps) this action will allow the merge of these developments with
  other developments of the year, building the future NEMO.
