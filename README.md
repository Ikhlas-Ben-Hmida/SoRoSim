SoRoSim v 2.32

This version of the toolbox allows the modelling and simulation of static and dynamic soft, open-chain robots.

You can start using this toolbox by defining links in MATLAB's command window, in order to create a link
define a variable as shown below:

LinkName = Link

You will then have to answer questions in the form of dialougue boxes in order to define the link's gemetric and material properties. After your link is created you can create a linkage which connects any links you've created previously. A linkage is a chain of 1 or more links 
connected in series to form an open chain and can be defied as shown below:

LinkageName = Linkage(LinkName1,LinkName2.....LinkNameN)

Dialougue boxes will then be used to collect input that is used to define
the linkage properties. You can access any link or linkage properties by calling the property name as shown below:

LinkName.LinkPropertyName

LinkageName.LinkagePropertyName

Similarly methods can also be used to evaluate properties outside of the classconstructor after the object is created 

LinkageName.MethodName(Input1...InputN)

You can run a static simulation by using the below command:
 
LinkageName.statics

the static simulation's output is a vector of joint angles
A dynamic simulation can be performed using:

LinkageName.dynamics

The two outputs of a dynamic simulation are the time vector (t) and a matrix containing the configuration (q) and velocity (qd) at every time element of (t). It is a good practice to save the outputs of the dynamic simulation for easy access.

The examples folder of the toolbox contains some saved linkages and links you can run simulations for.

More details about the theorey behind the toolbox and some of its applications can be found on:
https://arxiv.org/abs/2107.05494
