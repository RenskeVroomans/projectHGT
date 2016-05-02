# projectHGT

Recombination and horizontal gene transfer in the light of quasispecies theory, or something like that.
In this project, we can include the grind files, the python code, and potential individual-based models for easy collaboration.
Rejoice, we no longer have to copy files back and forth!

Bacio!

Important note: in order for this to work properly, you need a local git repository. In that way, you can work locally on your code,
and run it, and still have github keep track of our join progress.

How to do this:

make a local git repository:
	mkdir local_project
	cd local_project
	git init
	git config --global user.name "my name"
	git config --global user.email "my.address@gmail.com"


How to import the Github repository:
	git pull https://github.com/MyUserName/projectHGT
Add the github repository as a remote origin:
git remote add origin https://github.com/MyUserName/projectHGT

Now you can make changes to the files!
locally, you can keep track of your changes by committing them to your local repository:
	git commit -a -m "a handy message of the changes you made"

once you are ready to commit them to a branch or the master on github:
git push origin branch
