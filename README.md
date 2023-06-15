# ReproTox Paper Code

The Makefile has recipes for fetching all requisite data, constructing intermediary tables, and assembling all figures used in the paper, it can be viewed manually or used with the `make` command.

A makefile is of the form:
```
{thing_to_make}: {what is needed to make it}
  {commands to make the thing}
```
It can be used like `make {some_thing}` and it will be made along with any prerequisites.
