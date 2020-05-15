from HTMLParser import HTMLParser
import sys

# create a subclass and override the handler methods
class MyHTMLParser(HTMLParser):
    record = 0
    def handle_starttag(self, tag, attrs):
        if tag == "a":
            self.record += 1
    def handle_endtag(self, tag):
        if tag == "a":
            self.record -= 1
    def handle_data(self, data):
        if self.record > 0:
            self.value = data
