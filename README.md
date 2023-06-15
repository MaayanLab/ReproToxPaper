# ReproTox Paper Code

The Makefile has recipes for fetching all requisite data, constructing intermediary tables, and assembling all figures used in the paper, it can be viewed manually or used with the `make` command.

A makefile is of the form:
```
{thing_to_make}: {what is needed to make it}
  {commands to make the thing}
```
It can be used like `make {some_thing}` and it will be made along with any prerequisites.

## Figure 1
Figure 1 is a manually created graphic, as such there is no code.

## Figure 2-4
```bash
# see Makefile for details
make fig-2
make fig-3
make fig-4
```

## Figure 5
Figure 5 is a screenshot from the [ReproToxKG Website](https://maayanlab.cloud/reprotox-kg).