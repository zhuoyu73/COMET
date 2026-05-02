# MRFIL Reconstruction Library

Contains core routines to support iterative MR reconstruction in tandem with the Iterative Reconstruction Toolbox (IRT) maintained by Jeff Fessler at the University of Michigan.

## Installation

Pull a copy with your favorite git client

## Usage

To configure MATLAB to use our library, use the following lines in your startup.m file to set your paths correctly:

```MATLAB

% Replace [PATH_TO] with the path where the directories and repositories can be found.
basepath ='[PATH_TO]';

%Script to add necessary paths to repository
initializePaths;

```


## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a merge request to distribute the changes in your branch to the rest of the group. :-D

## History
2017-Sep-12 AMC: Beating back IRT with a weed whacker to critical routines only. These files will diverge with the distributed version of the IRT.
2016-Jul-07 AMC: Trimming repository to focus only on core MRFIL code and eliminate redundancies with the IRT.

## Credits

TODO: Write credits
