

# Keywords between dollar signs are substituted by subversion.
# The date and build will only tell the truth after a branch or tag,
# since different files in trunk will have been changed at different times
date ="$Date: 2012-01-31 12:19:17 -0800 (Tue, 31 Jan 2012) $".split()[1]
revision = "$Revision: 132 $".split()[1]


__version__ = '3.3' 


description = "CoreBio %s (%s)" % (__version__,  date)


