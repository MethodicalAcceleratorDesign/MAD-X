#/usr/bin/python
import smtplib
from email.MIMEText import MIMEText

def notify(recipient,subject,message):
    # Set up a MIMEText object (it's a dictionary)

    if recipient == 'jean-luc':
        toAdrs = 'Jean-Luc.Nougaret@cern.ch'
    elif recipient == 'admin':
        toAdrs = 'mad-automation-admin@cern.ch'
    else:
        toAdrs = 'Jean-Luc.Nougaret@cern.ch'
    
    TEXT = message

    msg = MIMEText(TEXT)
    
    msg['Subject'] = subject   
    msg['To'] = toAdrs
    # if above missing, message still gets delivered
    # but with 'undisclosed recipients'

    FROM = "Jean-Luc.Nougaret@cern.ch"
    TO = [ toAdrs ] # must be a list

    SERVER = 'localhost'
    
    server = smtplib.SMTP(SERVER)
    server.sendmail(FROM,TO,msg.as_string())
    server.quit()
