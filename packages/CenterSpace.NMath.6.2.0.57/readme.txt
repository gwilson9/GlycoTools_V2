RELEASE NOTES 
CenterSpace.NMath 6.2.0
December 2019
CenterSpace Software, LLC

=======================================================
Contents
-------------------------------------------------------
1.  Introduction
2.  CenterSpace Contact Information
3.  Dependencies
4.  Namespaces
5.  Documentation
6.  Code Examples
7.  Licensing
8.  Technical Support
9.  Consulting Services
10. Upgrades


=======================================================
1.  Introduction
-------------------------------------------------------
Package CenterSpace.NMath contains foundational classes
for object-oriented numerics on the .NET platform.
CenterSpace.NMath includes complex number classes,
general vector and matrix classes, structured sparse
matrix classes and factorizations, general sparse matrix
classes and factorizations, general matrix decompositions,
least squares solutions, random number generators,
Fast Fourier Transforms (FFTs), numerical integration
and differentiation methods, function minimization,
curve fitting, root-finding, linear and nonlinear
programming.


=======================================================
2.  CenterSpace Contact Information
-------------------------------------------------------
CenterSpace Software, LLC
622 NW 32nd Street
Corvallis, OR 97330 
USA
+1.541.896.1301

Web:          http://www.centerspace.net
Information:  info@centerspace.net
Sales:        sales@centerspace.net
Support:      support@centerspace.net
Skype:        centerspace


=======================================================
3. Dependencies
-------------------------------------------------------

o  Microsoft.Solver.Foundation
   https://www.nuget.org/packages/Microsoft.Solver.Foundation
   

=======================================================
4. Namespaces
-------------------------------------------------------
CenterSpace.NMath types are organized into three namespaces:

o  CenterSpace.NMath.Core

   Contains basic vector and matrix classes, complex number
   classes, and random number generators.
   
o  CenterSpace.NMath.Matrix

   Contains structured sparse matrix classes and factorizations,
   general sparse matrix classes and factorizations, general
   matrix decompositions, least squares solutions, and solutions
   to eigenvalue problems.
   
o  CenterSpace.NMath.Analysis

   Contains optimization, root-finding, and linear programming.
   
To avoid using fully qualified names, preface your code with the
appropriate namespace statements.


=======================================================
5.  Documentation
-------------------------------------------------------

o  User’s Guide
   http://www.centerspace.net/doc/NMath/NMath.UsersGuide.pdf
   http://www.centerspace.net/doc/NMath/user/
   
o  Reference Guide
   http://www.centerspace.net/doc/NMathSuite/ref/html/N_CenterSpace_NMath_Analysis.htm
   
o  Benchmarks
   http://www.centerspace.net/doc/NMath/whitepapers/NMath.Benchmarks.pdf
   
o  Whitepaper
   http://www.centerspace.net/doc/NMath/whitepapers/NMath.Whitepaper.pdf
   
o  Changelog
   http://www.centerspace.net/doc/NMath/changelog.txt


=======================================================
6.  Code Examples
-------------------------------------------------------
Numerous code examples are available for CenterSpace.NMath.

http://www.centerspace.net/examples/nmath/

Studying the examples is one of the best ways to learn
how to use the package.


=======================================================
7.  Licensing
-------------------------------------------------------
CenterSpace.NMath is distributed as a free, fully-functional
trial version for a 30-day evaluation period. For continued
use, please purchase a license through our secure online
store.

http://www.centerspace.net/order

CenterSpace.NMath is licensed per developer seat. Each
developer writing code using CenterSpace.NMath requires
a license. (For example, if 5 developers are writing an
application, and only 3 are doing most of the intense math
work, all 5 must have a license.) There are no runtime or
deployment fees for products developed that make use of
CenterSpace packages. The complete CenterSpace.NMath
license agreement is available here.

http://www.centerspace.net/license-agreement

CenterSpace.NMath license information is stored in a
license key which must be found at runtime. The license
key governs the properties of your installation. If no
license key is found, a default evaluation license key
is used which provides a free 30-day evaluation period
for CenterSpace.NMath on the current machine.

When you purchase one or more developer seats of
CenterSpace.NMath, you will be issued a license key
describing the terms of your license. To enter your
license key:

1. Open a NuGet console, by selecting Tools | NuGet
   Package Manager | Package Manager Console from 
   Visual Studio.
2. Type: NMathLicensing.exe
3. When the licensing tool launches, enter your name,
   email, and license key, and click OK.

The license key will be written to the registry.

You can also specify your license key using various other
mechanisms: by environment variable, by configuration app
setting, and programmatically. These mechanisms may be
preferable in group development environments, and at
deployment. See here for more information:

http://www.centerspace.net/doc/NMath/user/overview-83549.htm


=======================================================
8.  Technical Support
-------------------------------------------------------
CenterSpace technical support is available to all customers
currently under a maintenance contract. Product purchase
comes with one year of maintenance. To obtain technical
support, contact CenterSpace by email at:

support@centerspace.net

You can save time if you isolate the problem to a small
test case before contacting Technical Support. You may
also find your questions answered in our FAQ.

http://www.centerspace.net/faq


=======================================================
9.  Consulting Services
-------------------------------------------------------
CenterSpace developers are experts at object-oriented
numerics. We can bring valuable experience to companies
that are building and deploying critical .NET numerical
applications. Our developers can assist you in all
phases of your numerical development projects,
providing expertise and hands-on support in design,
development, optimization, and deployment. For more
information, contact sales@centerspace.net.


=======================================================
10.  Upgrades
-------------------------------------------------------
Upgrades are available free of charge to customers with
current annual maintenance contracts.
