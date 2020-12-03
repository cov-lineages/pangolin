### Notes on Pangolin Testing

The pangolin repository uses [Github Actions](https://github.com/features/actions) to provide continuous
integration testing. The files for this are in [.github/workflows].

#### Testing philosophy

The [.github/workflows/pangolin.yml] workflow does an install of pangolin followed by a simple run to
test that the tool is working. Tests can be added to that file to validate the output of the pangolin command.
This is essentially a functional test for the tool.

The test currently runs on every single pull request or push. It could be optimised, for example, to only
respond to action on some branches, or to do a check of which files have been changed before firing off the
complete testing workflow. For example, one should not need to test the pangolin tool if the only thing that
has changed is the README.md.

Also left as work for the future is a unit testing framework to test the components that make up pangolin.
