# Contributing to this project
Welcome, and thanks for your intrest in contributing to the project. Below are a few guidelines to help you contribute in a productive way. These are guidelines, not strict rules. If in doubt, use common sense and feel free to suggests changes to these guidelines via a pull request.

## Reporting bugs
Something not working right? Log it on the issue tracker. We'll have a much better
chance of fixing it if you can provide a reproducible example and as much detail
about the problem as possible.

## Addressing open issues
The issue tracker contains a collection of issues that are waiting for someone to handle them. A number of
labels are used to categorise these:

|  Label           |    Description |
| -------          |  --------------- |
| `bug`            |  Something is broken and should be fixed sooner rather than later. |
| `enhancement`    |  A new feature or substantial improvement to the algorithm. |
| `infrastructure` |  Tasks related to the general infrastructure of the project. These don't require any changes to the underlying algorithm. |
| `investigate`    |  Issues that may require some digging into the data and the output the algorithm generates from it to decide whether a change to the algorithm is required. This may result in a request for a new feature, a tweak to the algorithm or a bug report, rather than any code being contributed directly. |
| `question`       | Requests for clarification or discussion that aren't feature requests or bug reports. |
| `tweak`          |  A minor change to the algorithm or its implementation. |

If the issue description is unclear don't hesitate to ask for clarification. Once
you've identified a problem that you'd like to tackle, leave a comment to say so to avoid duplication of effort.

Use [feature branches](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow)
to work on new features and bug fixes. Feel free to create a pull request for the
branch at anytime to initiate discussion. It is helpful to label pull requests for
branches that are not yet ready for merging as `work-in-progress`. When a pull request
is created on GitHub it will trigger automatic checks by Travis CI, Coveralls and Landscape
to ensure everything still works and is looking good. It is advisable to look at the output
of these and fix any issues that may exist. 

## Contributing new features
Got an idea for a new feature or an improvement for the core algorithm?
We want to hear about it! Before starting work on an implementation it is 
advisable to discuss your idea. Open an issue with a description of what you want
to do and solicit comments.

## Styleguides
### Code formatting
Use **four spaces** for indenting.

Follow the standard Python [style for code](https://www.python.org/dev/peps/pep-0008/).
The easiest way to do this is to run [pylint](https://www.pylint.org/). Plugins for
many popular text editors and IDEs are available.

### Documentation
Use [Google style Python docstrings](http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html)
to document modules, classes and functions.

### Producing messages during program execution

* Use Python's [logging facilities](https://docs.python.org/2/howto/logging.html) 
to print messages during program execution. 
* Avoid the use of module level loggers, just create one when needed.
* Use `__name__` as name for loggers, except inside functions used as command-
  line scripts where you should use the package name instead. This will automatically
  generate an appropriate hierarchy of loggers.

