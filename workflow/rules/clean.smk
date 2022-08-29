rule cleanAfterAlignment:
    shell: """
find results/*/*/*/  -maxdepth 0 -not -name "bam" -exec rm -rf -v {{}} \;
rm -rf -v results/trackhub    
"""
