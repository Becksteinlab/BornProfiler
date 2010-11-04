#!/usr/bin/python
"""
:Author: Kaihsu Tai
:Year: 2008
:Licence: GPL
:Copyright: (c) 2008 Kaihsu Tai
:URL: http://en.wikiversity.org/wiki/Talk:Poisson%E2%80%93Boltzmann_profile_for_an_ion_channel

I wrote some Python code to automate this process. The job submission requires a queuing system called Grid Engine. Copyright © 2008 Kaihsu Tai. Moral rights asserted (why?). Hereby licensed under either GFDL or GNU General Public License at your option.

Obviously, you will have to change the parameters in the driver-functions (main()) to fit the purposes of this tutorial. That is the exercise for the reader! -- Kaihsu 16:12, 23 December 2008 (UTC)

Also, you will need to generate the PQR file, along with a file containing a list of the coordinates for the sample points. The script submits the job to a queueing system called Grid Engine, but you can submit the job by hand if you do not have this installed. -- Kaihsu 11:04, 9 January 2009 (UTC)

Additionally, you will need to change the ion species. -- Kaihsu 17:03, 30 January 2009 (UTC)

See also Two-dimensional Poisson-Boltzmann profile for an ion channel. -- Kaihsu 13:07, 30 April 2009 (UTC)
"""
import Gnuplot, os
 
class analyze:
  "analyze APBS energy profiling results"
 
  def readPoints(self):
    pointsFile = open(self.pointsName, "r")
    lines = pointsFile.readlines()
    pointsFile.close()
 
    points = []
    for line in lines:
      tokens = line.split()
      parsed = [float(tokens[0]), float(tokens[1]), float(tokens[2])]
      points.append(parsed)
    self.points = points
 
  def accumulate(self):
    self.readPoints()
    self.zE = []
    for point in self.points:
      z = point[2]
      outName = self.jobName + "/job_" + str(z) + ".out"
      lines = ""
      if (os.path.exists(outName)):
        outFile = open(outName, "r")
        lines = outFile.readlines()
        outFile.close()
      for line in lines:
        if (line[0:18] == "  Local net energy"):
          self.zE.append([z, float(line.split()[6])])
 
  def write(self):
    outName = self.jobName + ".dat"
    outFile = open(outName, "w")
    outFile.write("# z/angstrom E/(kJ/mol)\n")
    for each in self.zE:
      outFile.write("%8.3f %8.3e\n" % (each[0], each[1]))
    outFile.close()
 
  def plot(self):
    outName = self.jobName + ".dat"
    plotName = self.jobName + ".eps"
 
    # initialize
    g = Gnuplot.GnuplotProcess()
    cmd  = "set terminal postscript eps colour\n"
    cmd += 'set output "' + plotName + '"\n'
    cmd += """set style data lines
set xlabel "z / nm"
set ylabel "energy / (kJ/mol)"
"""
    cmd += 'plot "' + outName + '" using ($1/10):($2) title "' + self.jobName + '"\n'
 
    # do it
    g(cmd)
 
  def run(self):
    self.accumulate()
    self.write()
    self.plot()
 
def main():
  myP = analyze()
  myP.pointsName = "/sansom/s66/kaihsu/works/oxford/physiome/magnesium/samplepoints011.dat"
  myP.jobName = "H011cytoBut"
  myP.run()
 
if __name__ == "__main__":
  main()
